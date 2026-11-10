# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Run BioSQL tests using MySQL through the pymysql driver."""

import unittest

# Really do want "import *" to get all the test classes:
from common_BioSQL import *  # noqa: F403
from common_BioSQL import check_config

# Import these explicitly to avoid flake8 F405 below:
from common_BioSQL import load_biosql_ini

from BioSQL import BioSeqDatabase
from BioSQL.BioSeqDatabase import MysqlConnectorAdaptor
from BioSQL.DBUtils import Mysql_dbutils

DBDRIVER = "pymysql"
DBTYPE = "mysql"

DBHOST, DBUSER, DBPASSWD, TESTDB = load_biosql_ini(DBTYPE)

# This will abort if driver not installed etc:
check_config(DBDRIVER, DBTYPE, DBHOST, DBUSER, DBPASSWD, TESTDB)


class PymysqlConnectionTest(unittest.TestCase):
    """Check pymysql gets the MySQL-specific handling.

    The test classes imported from common_BioSQL read the driver from that
    module's globals, which the last BioSQL test module imported overwrites.
    This class uses this module's own settings, so it always talks to MySQL
    through pymysql.
    """

    def setUp(self):
        self.server = BioSeqDatabase.open_database(
            driver=DBDRIVER, host=DBHOST, user=DBUSER, passwd=DBPASSWD
        )

    def tearDown(self):
        self.server.close()

    def test_adaptor_and_dbutils(self):
        """Use the MySQL adaptor and last-insert-id helper."""
        self.assertIsInstance(self.server.adaptor, MysqlConnectorAdaptor)
        self.assertIsInstance(self.server.adaptor.dbutils, Mysql_dbutils)

    def test_ansi_quotes(self):
        """Set the ANSI_QUOTES SQL mode, so "rank" is a column name."""
        sql_mode = self.server.adaptor.execute_and_fetchall(
            "SELECT @@SESSION.sql_mode;"
        )[0][0]
        self.assertIn("ANSI_QUOTES", sql_mode.split(","))


if __name__ == "__main__":
    # Run the test cases
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
