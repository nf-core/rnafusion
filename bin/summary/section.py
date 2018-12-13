class Section:
    def __init__(self):
        self.id = ''
        self.title = ''
        self.subtitle = ''
        self.content = ''
        self.data = []
        self.graphs = {}
    
    def set_id(self, id):
        if id.strip():
            self.id = id.strip()

    def set_title(self, title):
        if title.strip():
            self.title = title.strip()
            
    def set_subtitle(self, subtitle):
        if subtitle.strip():
            self.subtitle = subtitle.strip()
            
    def set_content(self, content):
        if content.strip():
            self.content = content.strip()
    
    def set_query(self, db, query, query_param):
        self.data = db.select(query, query_param)

    def add_graph(self, graph):
        if graph.id not in self.graphs.keys():
            self.graphs[graph.id] = graph