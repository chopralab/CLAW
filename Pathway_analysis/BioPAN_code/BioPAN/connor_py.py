import psycopg2

# Connect to your postgres DB
conn = psycopg2.connect(
    dbname="biopan",
    user="biopan",
    password="access_biopan",
    host="localhost"
)

# Open a cursor to perform database operations
cur = conn.cursor()

# Read the SQL script
with open('database/create_biopan_tables.sql', 'r') as f:
    sql_script = f.read()

# Execute the SQL script
cur.execute(sql_script)

# Close communication with the database
cur.close()
conn.close()
