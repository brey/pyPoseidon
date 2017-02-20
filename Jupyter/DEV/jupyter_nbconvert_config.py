c = get_config()
#c.Exporter.preprocessors = [ 'bibpreprocessor.BibTexPreprocessor', 'pymdpreprocessor.PyMarkdownPreprocessor' ]
c.Exporter.preprocessors = [ 'pre_markdown.py.MarkdownPreprocessor' ]
c.Exporter.template_file = 'revtex_nocode.tplx'
