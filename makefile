SHELL = /bin/bash

install:
	@./check_moog.sh
	@tar -zxvf FASMA/models/apogee_kurucz.tar.gz
	@tar -zxvf FASMA/models/marcs.tar.gz
	@cp -r apogee_kurucz FASMA/models
	@cp -r marcs FASMA/models
	@rm -r apogee_kurucz
	@rm -r marcs
	@echo "Atmosphere models installed in dir: models"
	@echo "Installing dependencies..."
	@pip install -r requirements.txt
	@pip install .
	@echo "Dependencies installed"
	@echo ""
	@echo "FASMA is successfully installed!"
