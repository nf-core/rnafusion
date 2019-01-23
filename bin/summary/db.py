import sqlite3

class Db:
    def __init__(self, db_file=None):
        self.connect(db_file)

    def connect(self, db_file):
        try:
            if db_file is None:
                db_file = '/script-db/fusions.db'
            self.connection = sqlite3.connect(db_file)
            self.connection.row_factory = self.__dict_factory
        except sqlite3.Error as error:
            exit(error)

    def select(self, query, query_params=None):
        try:
            cur = self.connection.cursor()
            if query_params is None:
                cur.execute(query)
            else:
                cur.execute(query, query_params)
            res = cur.fetchall()
            cur.close()
            return res
        except sqlite3.Error as error:
            exit(error)

    @classmethod
    def __dict_factory(cls, cursor, row):
        tmp_dictionary = {}
        for idx, col in enumerate(cursor.description):
            tmp_dictionary[col[0]] = row[idx]
        return tmp_dictionary

    def __exit__(self, *args):
        self.connection.close()
