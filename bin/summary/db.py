import sqlite3

class Db:
    def __init__(self, db = None):
        self.connect(db)
    
    def connect(self, db):
        try:
            if db is None:
                db = '/script-db/fusions.db'
            self.connection = sqlite3.connect(db)
            self.connection.row_factory = self.__dict_factory
        except sqlite3.Error as error:
            exit(error)

    def select(self, query, query_params = []):
        try:
            cur = self.connection.cursor()
            if query_params is []:
                cur.execute(query)
            else:
                cur.execute(query, query_params)
            res = cur.fetchall()
            cur.close()
            return res
        except sqlite3.Error as error:
            exit(error)

    def __dict_factory(self, cursor, row):
        d = {}
        for idx, col in enumerate(cursor.description):
            d[col[0]] = row[idx]
        return d

    def __exit__(self, *args):
        self.connection.close()
