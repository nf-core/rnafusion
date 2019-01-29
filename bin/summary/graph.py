"""Module for Graph"""
class Graph:
    """Interface for Graph"""
    def __init__(self, graph_id, title, subtitle, data):
        self.graph_id = graph_id
        self.title = title
        self.subtitle = subtitle
        self.graph = data
