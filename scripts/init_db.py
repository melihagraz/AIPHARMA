from sqlalchemy import inspect
from app.db import engine, Base
import app.models  # modelleri yükle ki tablolar bilinsin

Base.metadata.create_all(bind=engine)

insp = inspect(engine)
print("Tables:", insp.get_table_names())
