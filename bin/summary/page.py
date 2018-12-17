from .master_page import MasterPage

class Page(MasterPage):
    def __init__(self, title, page_variables, partial_template):        
        self.page_variables = page_variables
        self.sections = {}
        super().__init__(title, partial_template)
    
    def add_section(self, section):
        if section.id not in self.sections.keys():
            self.sections[section.id] = section
        
    def get_section(self, section_id):
        return {} if self.sections[section_id] is None else self.sections[section_id]

    def get_content(self):
        return vars(self)