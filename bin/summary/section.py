class Section:
    def __init__(self):
        self.section_id = ''
        self.title = ''
        self.subtitle = ''
        self.content = ''
        self.data = []
        self.graphs = {}

    def set_id(self, section_id):
        if section_id.strip():
            self.section_id = section_id.strip()

    def set_title(self, title):
        if title.strip():
            self.title = title.strip()

    def set_subtitle(self, subtitle):
        if subtitle.strip():
            self.subtitle = subtitle.strip()

    def set_content(self, content):
        if content.strip():
            self.content = content.strip()

    def set_data(self, data):
        if data.strip():
            self.data = data

    def add_graph(self, graph):
        if graph.graph_id not in self.graphs.keys():
            self.graphs[graph.graph_id] = graph
