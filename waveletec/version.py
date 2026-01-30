__version__ = "0.3.0.1.0+nee_et.0.1.0"
full_version = __version__

release = 'dev' not in __version__ and '+' not in __version__
short_version = __version__.split("+")[0]
