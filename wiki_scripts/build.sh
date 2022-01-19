python3.8 -m venv venv38
pip install -r requirements.txt
source venv38/bin/activate
mkdocs build
htmlark site/print_page/index.html -o standalone.html
