# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Run BioSQL tests using MySQL through the pymysql driver."""

import unittest

import support

# Really do want "import *" to get all the test classes:
from common_BioSQL import *  # noqa: F403
from common_BioSQL import check_config

# Import these explicitly to avoid flake8 F405 below:
from common_BioSQL import load_biosql_ini

from Bio import SeqIO
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


class PymysqlLoadTest(unittest.TestCase):
    """Load the schema and records through pymysql, on this module's settings.

    Uses a database of its own, so it cannot collide with the shared tests of
    the other MySQL modules on the same server.
    """

    testdb = TESTDB + "_pymysql"

    def _drop_database(self):
        server = BioSeqDatabase.open_database(
            driver=DBDRIVER, host=DBHOST, user=DBUSER, passwd=DBPASSWD
        )
        try:
            server.adaptor.execute(f"DROP DATABASE IF EXISTS {self.testdb}")
        finally:
            server.close()

    def setUp(self):
        self._drop_database()
        server = BioSeqDatabase.open_database(
            driver=DBDRIVER, host=DBHOST, user=DBUSER, passwd=DBPASSWD
        )
        try:
            server.adaptor.execute(f"CREATE DATABASE {self.testdb}")
        finally:
            server.close()
        self.addCleanup(self._drop_database)

    def test_load_schema_and_records(self):
        """Load the BioSQL schema, then GenBank records with ranked features."""
        server = BioSeqDatabase.open_database(
            driver=DBDRIVER,
            host=DBHOST,
            user=DBUSER,
            passwd=DBPASSWD,
            db=self.testdb,
        )
        self.addCleanup(server.close)
        # Raised "Module pymysql not supported by the loader." before:
        server.load_database_sql(support.DATA / "BioSQL" / "biosqldb-mysql.sql")
        db = server.new_database("pymysql_load")
        records = list(SeqIO.parse(support.DATA / "GenBank" / "cor6_6.gb", "gb"))
        # Without ANSI_QUOTES, MySQL rejects the quoted "rank" column (1064):
        self.assertEqual(db.load(records), len(records))
        server.commit()
        for record in records:
            stored = db.lookup(accession=record.id.split(".")[0])
            self.assertEqual(stored.seq, record.seq)
            self.assertEqual(len(stored.features), len(record.features))


if __name__ == "__main__":
    # Run the test cases
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
