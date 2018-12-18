import os
from pathlib import Path
from .report_config import ReportConfig
from jinja2 import Environment, FileSystemLoader, Markup

class Report:
    def __init__(self, config, output_dir):
        self.j2_env = Environment(
            loader=FileSystemLoader(os.path.dirname(os.path.abspath(__file__))),
            trim_blocks=True,
            autoescape=True
        )
        self.pages = {}
        self.output_dir = output_dir
        self.config = ReportConfig(config)
        self.j2_variables = self.config.get_variables()

        if not os.path.exists(self.output_dir):
            os.mkdir(output_dir)

        # Helper fusions
        self.j2_env.globals['include_raw'] = self.__include_raw

    def add_page(self, page):
        if page.filename in self.pages:
            print('Page ' + page.filename + ' already exists, skipping ...')
        else:
            dict_page = page.get_content()
            template_variables = {**self.j2_variables, **dict_page}
            template_variables['menu'] = [(key, page.get_section(key).title) for key in page.sections.keys()]
            self.pages[page.filename] = template_variables

    def render(self):
        for page, template_variables in self.pages.items():
            # print('Building ' + page +  ' ...')
            output = self.j2_env.get_template('template.html').render(template_variables)
            with open(os.path.join(self.output_dir, page), 'w') as file_out:
                file_out.write(output)

    def __include_raw(self, filename):
        file_extension = Path(filename).suffix
        if file_extension == '.css':
            return Markup('<style type="text/css">' + self.j2_env.loader.get_source(self.j2_env, filename)[0] + '</style>')
        if file_extension == '.js':
            return Markup('<script>' + self.j2_env.loader.get_source(self.j2_env, filename)[0] + '</script>')

        return Markup(self.j2_env.loader.get_source(self.j2_env, filename)[0])