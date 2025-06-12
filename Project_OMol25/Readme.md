1) You will require Postgres Database, of 4M data. (first requirement)
2) In that database, the column ("chemical coordinates") have been hashed. We are using that hash column only, for the search query based on 4M chemical componds data. 
3) You can see the configuration for the database, in the end of the server.py
4) Once the setup is completed. 
5) Create a python enviroment, and enter that python enviroment. (commands differ in windows and linux)
6) In the same terminal go to ```Project_OMol25/server``` and run ```python bridge.py```. 
7) Now open another terminal and go to directory ```Project_OMol25/client``` in the terminal. 
8) Run ```npm i```, ```npm run dev``` in the terminal.
9) Your project must be running on ```localhost:5173```