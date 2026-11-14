# AIPHARMA araştırma yazılımı yayın adayı — kullanıcı adımları

Bu ZIP, mevcut `~/AIPHARMA` projesine **eklenen** dosyalardan oluşur; proje kodunu veya veri sonuçlarını silmez. Bu paket tek başına tam AIPHARMA kaynak kodu değildir. Projenin kendi WSL kopyasından güvenli, sürümlü bir source candidate üretir.

```bash
cd ~/AIPHARMA
source .venv/bin/activate
unzip -n /mnt/c/Users/cobre/Downloads/aipharma_releasekit_rc1.zip -d ~/AIPHARMA
python scripts/aipharma_release.py inspect
python scripts/aipharma_release.py freeze
python scripts/aipharma_release.py test
python scripts/aipharma_release.py bundle
```

Oluşacak dosyalar:
- `dist_release/AIPHARMA_0.1.0-rc1_source_candidate.zip`
- `dist_release/AIPHARMA_0.1.0-rc1_source_candidate.zip.sha256`
- `dist_release/release_audit.json`
- `reproducibility/local/environment.json` (yerel, kamuya açık değildir)
- `reproducibility/local/pip-freeze-pyXY.txt` (kontrol edilmeli)

`inspect` çıktısında public repository, LICENSE veya model/data eksiklikleri `publication_blockers` olarak görünmesi beklenir. Bunları otomatik olarak var saymıyoruz. GitHub/Zenodo'ya yükleme yapılmaz.

**Dikkat:** `bundle` kaynak kodu küçük, önceden belirlenmiş uzantı ve klasörlerden toplar, bilimsel veri ve ağırlıkları dışarıda bırakır. Yayınlamadan önce arşivi manuel olarak açıp fikrî mülkiyet ve güvenlik kontrolü yap. Kök lisans için hak sahiplerinden izin al.

**Eski README uyarısı:** Önceki kök `README.md` Phase 6 ifadeleri içerebilir. `reproducibility/README_RELEASE_DRAFT.md` güncel sonuçları anlatan bir taslaktır. İçeriği gözden geçirip kamuya açmadan önce kök README'ye taşı; otomatik olarak orijinal README'nin üzerine yazmıyoruz.
