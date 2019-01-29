"""Module for Section"""
class Section:
    """This class is used to define a Section in a report with predefined variables"""
    def __init__(self):
        self.section_id = ''
        self.title = ''
        self.subtitle = ''
        self.content = ''
        self.data = []
        self.graphs = {}

    def set_id(self, section_id):
        """Helper method for setting section id"""
        if section_id.strip():
            self.section_id = section_id.strip()

    def set_title(self, title):
        """Helper method for setting section title"""
        if title.strip():
            self.title = title.strip()

    def set_subtitle(self, subtitle):
        """Helper method for setting section subtitle"""
        if subtitle.strip():
            self.subtitle = subtitle.strip()

    def set_content(self, content):
        """Helper method for setting section content"""
        if content.strip():
            self.content = content.strip()

    def set_data(self, data):
        """Helper method for setting data content, usually generated from local database but can
        be custom"""
        if data.strip():
            self.data = data

    def add_graph(self, graph):
        """Helper method for adding graph to the section, unlimited number"""
        if graph.graph_id not in self.graphs.keys():
            self.graphs[graph.graph_id] = graph
