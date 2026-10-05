"""scripts/build_hpo_index.py: the hp.index builder that scripts/fetch_data.sh calls. No network, no FastHPOCR run."""
import pytest
from scripts import build_hpo_index


def test_missing_obo_exits_with_download_hint(tmp_path):
    with pytest.raises(SystemExit) as e:
        build_hpo_index.main(["--obo", str(tmp_path / "nope.obo"), "--out-dir", str(tmp_path)])
    msg = str(e.value)
    assert "hp.obo not found" in msg and "purl.obolibrary.org/obo/hp.obo" in msg
