"""Module for Report"""
import os
from pathlib import Path
from jinja2 import Environment, FileSystemLoader, Markup
from .report_config import ReportConfig

class Report:
    """Class containing core functionality for generating report"""
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

        # Making sure output directory exists
        if not os.path.exists(self.output_dir):
            os.mkdir(output_dir)

        # Helper fusion for including raw content in Jinja
        self.j2_env.globals['include_raw'] = self.__include_raw

    def add_page(self, page):
        """
        Helper function for adding a page
        Args:
            page (class Page)
        """
        if page.filename in self.pages:
            print('Page ' + page.filename + ' already exists, skipping ...')
        else:
            dict_page = page.get_content()
            template_variables = {**self.j2_variables, **dict_page}
            template_variables['menu'] = [
                (key, page.get_section(key).title) for key in page.sections.keys()
            ]
            # using filenames as identifiers
            self.pages[page.filename] = template_variables
            # after adding immediately generate fusion page
            self.render(page.filename, template_variables)

    def render(self, filename, template_variables):
        """Helper method rendering page using Jinja template engine"""
        output = self.j2_env.get_template('template.html').render(template_variables)
        with open(os.path.join(self.output_dir, filename), 'w') as file_out:
            file_out.write(output)

    def __include_raw(self, filename):
        """Helper fusion for including raw content in Jinja, mostly used to include custom
        or vendor javascript and custom css"""
        file_extension = Path(filename).suffix
        if file_extension == '.css':
            return Markup(
                '<style type="text/css">{css}</style>'.format(
                    css=self.j2_env.loader.get_source(self.j2_env, filename)[0]
                )
            )
        if file_extension == '.js':
            return Markup(
                '<script>{js}</script>'.format(
                    js=self.j2_env.loader.get_source(self.j2_env, filename)[0]
                )
            )

        return Markup(self.j2_env.loader.get_source(self.j2_env, filename)[0])
