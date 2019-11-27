SHELL = /bin/bash

install:
	@./check_moog.sh
	@mkdir -p FASMA/linelist
	@mkdir -p FASMA/spectra
	@tar -zxvf FASMA/models/apogee_kurucz.tar.gz
	@mv apogee_kurucz FASMA/models
	@tar -zxvf FASMA/models/marcs.tar.gz
	@mv marcs FASMA/models
	@echo "Atmosphere models installed in dir: models"
	@echo "Installing dependencies..."
	@pip install -r requirements.txt
	# @conda install -c anaconda wxpython=3.0.0.0
	@echo "Dependencies installed"
	@echo ""
	@echo "FASMA is successfully installed!"
	@echo "Type"
	@echo "    python FASMA.py"
	@echo "to start. Happy spectroscopying :)"
