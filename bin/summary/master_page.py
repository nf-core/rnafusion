class MasterPage:
    def __init__(self, title, partial_template):
        self.title = title
        self.filename = title.replace('--', '_') + ".html"
        self.dynamic_partial = 'partials/' + partial_template + ".html"

    def get_content(self):
        return vars(self)
