"""Module for Report Configuration"""
import os
import base64
from datetime import datetime
from yaml import safe_load, YAMLError


class ReportConfig:
    """Class for adjusting report"""
    def __init__(self, config=None):
        self.current_path = os.path.dirname(os.path.abspath(__file__))
        self.report_title = 'nfcore/rnafusion summary report'
        self.logo = base64.b64encode(
            open(os.path.join(self.current_path, 'assets/img/rnafusion_logo.png'), 'rb').read()
        ).decode('utf-8')
        self.institution = ''
        self.date_format = '%d/%m/%Y'
        self.date = datetime.now().strftime(self.date_format)
        self.assets = {}

        if config is not None:
            self.__parse(config)

    def get_variables(self):
        """Helper method returning all variables. Used for Jinja templating"""
        return vars(self)

    def __parse(self, config):
        """
        Method for parsing the configuration file
        Args:
            config (string): path to configuration file
        """
        try:
            with open(config, 'r') as in_file:
                try:
                    data = safe_load(in_file)
                    self.__set_title(data)
                    self.__set_institution(data)
                    self.__set_date_format(data)
                    self.__set_assets(data)
                except YAMLError as error:
                    print(error)
        except IOError as error:
            print(error)
            print('Configuration file not found, skipping ...')

    def __set_title(self, config):
        """Helper function for setting a custom title"""
        if config['report_title'] is not None:
            self.report_title = config['report_title'].strip()

    def __set_institution(self, config):
        """Helper function for adding an institution logo"""
        if config['institution'] is not None and os.path.exists(config['institution']):
            self.institution = base64.b64encode(
                open(os.path.join(self.current_path, config['institution']), 'rb').read()
            ).decode('utf-8')

    def __set_date_format(self, config):
        """Helper function for setting a custom date format"""
        if config['date_format'] is not None:
            self.date_format = config['date_format']
            self.date = datetime.now().strftime(self.date_format)

    def __set_assets(self, config):
        """Helper function for adding custom CSS or Javascript"""
        if config['assets'] is not None:
            for key, value in config['assets'].items():
                if key in ('css', 'js'):
                    self.assets[key] = [x for x in value if os.path.exists(x)]
