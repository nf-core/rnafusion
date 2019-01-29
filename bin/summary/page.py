"""Module for Page"""
from .master_page import MasterPage

class Page(MasterPage):
    """This class is used to define a custom Page in a report, inherits some defaults from Master"""
    def __init__(self, title, page_variables, partial_template):
        self.page_variables = page_variables
        self.sections = {}
        super().__init__(title, partial_template)

    def add_section(self, section):
        """
        Method for adding new section if it doesn't exist
        Args:
            section (class Section)
        """
        if section.section_id not in self.sections.keys():
            self.sections[section.section_id] = section

    def get_section(self, section_id):
        """
        Method for returning section if it exists
        Args:
            section_id (string, defined by user, not generated in any way)
        """
        return {} if self.sections[section_id] is None else self.sections[section_id]

    def get_content(self):
        """Helper method returning all variables. Used for Jinja templating"""
        return vars(self)
