"""Module for MasterPage"""
class MasterPage:
    """Class defining main MasterPage. Variables are inherited for all created Pages (class Page)"""
    def __init__(self, title, partial_template):
        self.title = title
        self.filename = title.replace('--', '_') + ".html"
        self.dynamic_partial = 'partials/' + partial_template + ".html"

    def get_content(self):
        """Helper method returning all variables. Used for Jinja templating"""
        return vars(self)
