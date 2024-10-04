# don't run the "app.run()" line in Jupyter. 
# Instead, copy the text of  this block and save as "first_app.py" after uncommenting the last 2 lines and 
# run it from a terminal or prompt with this command: "python app_test.py"
# The BP.csv file must be in the same folder


import pandas as pd
from sqlalchemy import create_engine, text
from flask import Flask, request, jsonify
from flask_cors import CORS
import json

def get_chemicals(bp_value):
    chemicals = {}
    df = pd.read_csv("data/BP.csv")
    
    # Create SQLite engine
    engine = create_engine('sqlite://', echo=False)
    
    # Save DataFrame to SQL database
    df.to_sql('chemical', con=engine, if_exists='replace', index=False)
    
    # Define column names (in case you need to map them later)
    col_names = ['index', 'compound_number', 'name', 'BP_C', 'BP_K', 'SMILES', 'MW']
    
    try:
        # Establish a connection using context manager
        with engine.connect() as connection:
            
            # Use text() to run the SQL query
            query = text("SELECT * FROM chemical WHERE BP_C = :bp_value")
            result = connection.execute(query, {'bp_value': bp_value})
            
            # Fetch all rows matching the boiling point condition as mappings (dictionary-like rows)
            rows = result.mappings().all()  # Use .mappings() to access rows by column names
            
            if not rows:
                print("No results found")
                return json.dumps({})
            
            # Prepare the dictionary of chemicals
            for row in rows:
                # `row` is now a dictionary-like object, so you can access values by column name
                chemical = {col: row[col] for col in col_names if col in row}
                chemicals[row['name']] = chemical  # Use 'name' as the key for each chemical

    except Exception as e:
        print(f"An error occurred: {e}")
        return json.dumps({})

    # Convert dictionary to JSON for output
    chemicals_out = json.dumps(chemicals, separators=(',', ':'))
    
    return chemicals_out

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})

@app.route('/api/chemical/<name>', methods=['GET','POST'])

def api_get_users(name):
    return jsonify(get_chemicals(name))

if __name__ == ('__main__'):
    app.run()
