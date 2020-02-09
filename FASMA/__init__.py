from .synthDriver import synthMethod

if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        cfgfile = sys.argv[1]
    else:
        cfgfile = 'config.yml'
    driver = synthMethod(cfgfile=cfgfile, overwrite=None)
    driver.synthdriver()
