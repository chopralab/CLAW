import psycopg2

def execute_sql_file(filename, db_connection):
    with open(filename, 'r') as sql_file:
        sql_script = sql_file.read()

    cursor = db_connection.cursor()
    try:
        cursor.execute(sql_script)
        db_connection.commit()
        print(f"Executed file: {filename}")
    except (Exception, psycopg2.DatabaseError) as error:
        print(f"Error while executing file: {filename}, the error is: {error}")
    finally:
        cursor.close()

# Replace these with your actual details
DATABASE_NAME = "biopan"
USER = "biopan"
PASSWORD = "access_biopan"
HOST = "localhost"
PORT = "5432"

connection = psycopg2.connect(
    dbname=DATABASE_NAME, 
    user=USER, 
    password=PASSWORD, 
    host=HOST, 
    port=PORT
)

# Execute the SQL files
execute_sql_file('../database/create_biopan_tables.sql', connection)
execute_sql_file('../database/insert_gene_table.sql', connection)

connection.close()