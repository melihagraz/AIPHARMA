"""HTTP integration regression for login/profile and private v1 job resources."""
from app.core.config import get_settings
from app.models.auth_user import AuthUser


def _user(client, email):
    p='correct-horse-battery-phrase-12'
    r=client.post('/api/auth/register',json={'email':email,'password':p,'full_name':'Researcher'})
    assert r.status_code==200,r.text
    r=client.post('/api/auth/login',json={'email':email,'password':p})
    assert r.status_code==200,r.text
    return {'Authorization':'Bearer '+r.json()['access_token']}


def test_login_profile_and_jobs_are_scoped(client, test_db, monkeypatch):
    monkeypatch.setenv('SESSION_SIGNING_KEY','offline-integration-tests-only-secret-at-least-32')
    get_settings.cache_clear()
    try:
        assert client.get('/api/v1/jobs').status_code==401
        a=_user(client,'alice@example.org')
        b=_user(client,'bob@example.org')
        x=client.post('/api/v1/jobs',json={'receptor_name':'ALK5','tenant_id':'bob'},headers=a)
        assert x.status_code==200,x.text
        aid=x.json()['id']
        assert x.json()['status']=='queued' # no fabricated completion without worker
        assert x.json()['tenant_id']!='bob'
        assert client.get('/api/v1/jobs',headers=a).json()[0]['id']==aid
        assert client.get('/api/v1/jobs',headers=b).json()==[]
        assert client.get('/api/v1/jobs/'+aid,headers=b).status_code==404
        assert client.delete('/api/v1/jobs/'+aid,headers=b).status_code==404
        assert client.get('/api/v1/jobs/'+aid,headers=a).status_code==200
        assert client.get('/api/auth/profile',headers=a).json()['email']=='alice@example.org'
        assert client.get('/api/auth/profile',params={'email':'alice@example.org'}).status_code==401
        assert client.delete('/api/v1/jobs/'+aid,headers=a).status_code==200
        assert client.get('/api/v1/jobs',headers=a).json()==[]
        row=test_db.query(AuthUser).filter(AuthUser.email=='alice@example.org').first()
        assert row.password.startswith('$argon2id$')
    finally:
        get_settings.cache_clear()


def test_agent_chat_denies_anonymous(client):
    assert client.post('/api/agent/chat',json={'message':'hello'}).status_code in {401,422}
