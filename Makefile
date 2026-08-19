PYTHON ?= python3

.PHONY: setup pipeline dashboard

setup:
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -r requirements.txt

pipeline:
	$(PYTHON) pipeline.py

dashboard: pipeline
	$(PYTHON) -m streamlit run dashboard.py --server.address=0.0.0.0
