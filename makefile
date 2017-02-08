SHELL = /bin/bash

install:
	@./check_moog.sh
	@mkdir -p linelist
	@mkdir -p spectra
	@tar zxf models/apogee_kurucz.tar.gz
	@tar zxf models/marcs.tar.gz
	@echo "Atmosphere models installed in dir: models"
	@echo "Installing dependencies..."
	@pip install -r requirements.txt
	@conda install -c anaconda wxpython=3.0.0.0
	@echo "Dependencies installed"
	@echo ""
	@echo "MOOGme is successfully installed!"
	@echo "Type"
	@echo "    python FASMA.py"
	@echo "to start. Happy spectroscopying :)"
