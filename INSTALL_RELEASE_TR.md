# AIPHARMA v1.0.0 — GitHub yayın adayı

**Durum (14 Kasım 2026):** Kaynak kodunun seçilmiş sürümü GitHub\u0027daki **private** `release/v1.0.0-candidate` dalında. Henüz public GitHub Release veya AIPHARMA Zenodo DOI oluşturulmadı.

## İçerik ve sınırlamalar

- Kaynak kodu ve bilimsel denetim betikleri mevcuttur; indirilen model ağırlıkları, ham benchmark verileri, tam React arayüzü ve kişisel raporlar bu pakette yoktur.
- `backend`, `Reinvent_App`, `Admet_App`, `Docking_App` ve `Pocket_App` farklı bilimsel servisler ve bağımlılıklar kullanabilir. Temiz main API kurulumu bütün bilimsel modellerin çalıştığını göstermez.
- Önceki yerel testte 179 test geçti, 3 test deselected kaldı; bu bağımsız tam-yığın doğrulama değildir.

## Depoyu inceleme

    gh repo clone melihagraz/AIPHARMA
    cd AIPHARMA
    git switch release/v1.0.0-candidate

Kurulumun daha önce denetlenen kısmı için `SETUP.md`, kapsam ve sınırlamalar için `REPRODUCIBILITY.md`, kaynak hakları için `PUBLICATION_STATUS.md` ve `reproducibility/THIRD_PARTY_INVENTORY_DRAFT.md` belgelerini oku.

**Public yayın için:** Hak sahipliği ve üçüncü taraf lisans denetimi, nihai Git ağacının gizli bilgi/veri taraması, kaynak arşiv SHA-256 ve metadata incelemesini tamamla. Ardından sürüm etiketini ve GitHub Release\u0027i yayımla, gerçek Zenodo DOI\u0027sini `CITATION.cff` içine ekle. Onay verilmeden depo görünürlüğünü değiştirme.
