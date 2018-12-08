from .master_page import MasterPage

class Page(MasterPage):
    def __init__(self, db, title, partial_template):
        self.db = db
        self.content = {}
        super().__init__(title, partial_template)
    
    def add_custom_section(self, key, title, content):
        if key in self.content.keys():
            print('Key: ' + key + ' already exists, skipping ...')
        else:
            self.content[key] = {
                'title': title,
                'data': content
            }

    def add_section(self, key, title, query, query_params):
        if key in self.content.keys():
            print('Key: ' + key + ' already exists, skipping ...')
        else:
            self.content[key] = {
                'title': title,
                'data': self.db.select(query, query_params)
            }

    def get_section(self, key):
        if key not in self.content.keys():
            print('Key: ' + key + ' doesn\'t exists, skipping ...')
            return {}

        return self.content[key]

    def add_graph(self, key, graph_id, graph):
        if key not in self.content.keys():
            print('Key: ' + key + ' doesn\'t exists, skipping ...')
        else:
            if 'graphs' not in self.content[key].keys():
                self.content[key]['graphs'] = {}

            if graph_id in self.content[key]['graphs']:
                print('Key: ' + key + ' already exists, skipping ...')
            else:
                self.content[key]['graphs'][graph_id] = graph
