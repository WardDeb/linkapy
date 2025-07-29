from importlib.metadata import version as importlibversion

project = 'linkapy'
author = 'WardDeb'
version = importlibversion("linkapy")
release = importlibversion("linkapy")

extensions = ['autoapi.extension']
language = 'en'
master_doc = 'index'
pygments_style = 'sphinx'
source_suffix = '.rst'

html_theme = 'sphinx_rtd_theme'
autoapi_dirs = ['../python']
autodoc_mock_imports = ["linkapy.linkapy"]
