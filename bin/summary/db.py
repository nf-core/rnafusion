"""Module for accessing local database"""
import sqlite3

class Db:
    """Database wrapper around sqlite3 for summary report"""
    def __init__(self, db_file=None):
        self.connect(db_file)

    def connect(self, db_file):
        """
        Wrapper around default connect function
        Args:
            db_file: local database file *.db
        """
        try:
            if db_file is None:
                db_file = '/script-db/fusions.db'
            self.connection = sqlite3.connect(db_file)
            self.connection.row_factory = self.__dict_factory
        except sqlite3.Error as error:
            exit(error)

    def select(self, query, query_params=None):
        """
        Wrapper around default fetch function
        Args:
            query (string): SQL statement
            query_params (list): list of all parameters, SQL statement should be sanitized
        """
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
        """Helper class for converting SQL results into dictionary"""
        tmp_dictionary = {}
        for idx, col in enumerate(cursor.description):
            tmp_dictionary[col[0]] = row[idx]
        return tmp_dictionary

    def __exit__(self, *args):
        """Close connection on exit"""
        self.connection.close()
