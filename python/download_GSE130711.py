from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from bs4 import BeautifulSoup
import os

file_path = '/Users/megan/Projects/SnapHiC-D/ext/Lee_Astro_samples.txt'  # Replace with the actual file path
filenames = []
with open(file_path, 'r') as file:
    lines = file.readlines()
    for line in lines:
        file_name = line.strip()  # Remove leading/trailing whitespace and newline characters
        filenames.append(file_name)

file_path = '/Users/megan/Projects/SnapHiC-D/ext/Lee_MG_samples.txt'  # Replace with the actual file path
with open(file_path, 'r') as file:
    lines = file.readlines()
    for line in lines:
        file_name = line.strip()  # Remove leading/trailing whitespace and newline characters
        filenames.append(file_name)

file_names = [os.path.basename(filename) for filename in filenames]
print(file_names[0])
# Launch a web browser (you need to have the appropriate web driver installed)
driver = webdriver.Chrome()

# Navigate to the webpage containing the button
driver.get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130711")

# Find the button element by its ID
button_element = driver.find_element(By.ID, "customDl")

# Simulate clicking the button
button_element.click()

# Wait for the page to load or any necessary actions to complete
# You may need to adjust the sleep duration depending on the webpage's behavior
import time
time.sleep(5)
# Get the HTML content of the webpage
html_content = driver.page_source

# Parse the HTML content using BeautifulSoup
soup = BeautifulSoup(html_content, "html.parser")

# Find all input elements with a specific class
# input_elements = soup.find_all("input", class_="customDlGroup")

# Wait for the checkbox elements to be clickable
wait = WebDriverWait(driver, 10)
checkbox_elements = wait.until(
    EC.presence_of_all_elements_located((By.CLASS_NAME, "customDlGroup"))
)

# Iterate through the checkbox elements and check the checkbox if the filename matches
for checkbox_element in checkbox_elements:
    parent_element = checkbox_element.find_element(By.XPATH, "./..")
    text_filename = parent_element.text.strip()
    if text_filename in file_names:
        checkbox_element.click()

# Close the web browser
driver.quit()