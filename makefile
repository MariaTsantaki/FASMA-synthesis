SHELL = /bin/bash

install:
	@tar -zxvf FASMA/models/apogee_kurucz.tar.gz
	@tar -zxvf FASMA/models/marcs.tar.gz
	@cp -r apogee_kurucz FASMA/models
	@cp -r marcs FASMA/models
	@rm -r apogee_kurucz
	@rm -r marcs
	@echo "Atmosphere models installed in dir: models"
	@echo "Installing dependencies..."
	@sudo pip install -r requirements.txt
	@sudo pip install .
	@echo "Dependencies installed"
	@echo ""
	@./check_moog.sh
	@echo "FASMA is successfully installed!"
