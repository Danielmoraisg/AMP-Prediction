from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from bs4 import BeautifulSoup as bs
import pandas as pd
import xerox
import time
from Bio import SeqIO
import os
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
import random
import string
from tqdm import tqdm


def set_chrome_config(headless = True):
    options = webdriver.ChromeOptions()
    options.add_argument("no-sandbox")
    options.add_argument("start-maximized")
    options.add_argument("enable-automation")
    options.add_argument("--disable-gpu")
    options.add_argument("--disable-extensions")
    if headless == True:
        options.add_argument('--headless')
    options.add_argument("--window-size=1920,1080")
    options.add_argument('--ignore-certificate-errors')
    options.add_argument('--disable-dev-shm-usage')
    options.add_argument("--proxy-server='direct://'")
    options.add_argument("--proxy-bypass-list=*")
    options.add_argument("--lang=en_US")
    options.add_argument("--disable-infobars")
    options.add_argument("--disable-browser-side-navigation")
    return options
    
def random_string(length):
    # With combination of lower and upper case
    result_str = ''.join(random.choice(string.ascii_letters) for i in range(length))
    # print random string
    return result_str

def portreports(file, driver_path, url = 'https://www.portoreports.com/stm'):
    with open(file) as f:
        fasta = f.read()
    driver = webdriver.Chrome(driver_path,chrome_options = set_chrome_config())
    driver.get(url)
    inputElement = driver.find_element_by_xpath('//*[@id="post-77"]/div/form/textarea')
    #xerox.copy(fasta)
    driver.execute_script("arguments[0].value = arguments[1];", inputElement,fasta)
    #inputElement.send_keys(Keys.CONTROL+ "v")
    inputElement.submit()
    content = driver.page_source
    soup = bs(content,features="lxml")
    driver.close()
    df = pd.read_html(soup.prettify())[0]
    df['Sequence Number'] = df['Sequence Number'] - 2
    return df.iloc[2: , :].reset_index(drop = True)

#function to generate bathces of sequences
def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]    
        
#function to skip sequences in a range
def skip(file,upper_limit=9999999, lower_limit = 0, remove = True):
    wanted = []
    with open(file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            l =len(record.seq)
            if l<upper_limit and l>lower_limit:
                wanted.append(record)
        SeqIO.write(wanted, "temp.fasta", "fasta")
        file = 'temp.fasta'
    with open(file) as f:
        fasta = f.read()
    print(f' {len(wanted)} sequences were kept')
    if remove:
        os.remove("temp.fasta")
    return fasta

def dbaasp(file, driver_path,problems = False, n_try = 23, url = 'https://dbaasp.org/prediction/general', skip_100 = False, skip_4 = False):
    #skipping potentialy big or small sequences
    if skip_100 or skip_4:
        if skip_100 and skip_4:
            fasta = skip(file,upper_limit=100, lower_limit = 4, remove = False)
        elif skip_100:
            fasta = skip(file, upper_limit=100, remove = False)
        else:
            fasta = skip(file, lower_limit = 4, remove = False)
        file = 'temp.fasta'
    #generate tempfasta of 20 aas and passing it into the dbaasp website
    with open(file) as handle:
        records = SeqIO.parse(handle, "fasta")
        first = True
        if first:
            df0 = pd.DataFrame(columns = ['Seq. ID','Class']) #starting the DataFrame
            first = False
        problems = []
        for record in batch(list(records), n_try):
            SeqIO.write(record, "small_temp.fasta", "fasta")
            with open('small_temp.fasta') as f:
                small_fasta = f.read()
            os.remove("small_temp.fasta")
            driver = webdriver.Chrome(driver_path)
            driver.get(url)
            inputElement = driver.find_element_by_xpath('/html/body/main/div[2]/div/textarea')
            xerox.copy(small_fasta)
            inputElement.send_keys(Keys.CONTROL+ "v")
            driver.find_element_by_xpath('/html/body/main/div[2]/div/button').click()
            time.sleep(2)
            content = driver.page_source
            soup = bs(content,features="lxml")
            driver.close()
            try:
                df = pd.read_html(soup.prettify())[0]
                df.drop(df.tail(1).index,inplace=True)
                df0 = pd.concat([df0, df], ignore_index=True)
            except IndexError:
                problems.append(record)
        print(f'{len(problems)*23} were not able to be used by dbaasp\nTo check sequences that were not used please use the argument problems = true')
    os.remove('temp.fasta')
    if problems:
        return df0, problems
    else:
        return df0
    
def campr3(file, driver_path, url = 'http://www.camp.bicnirrh.res.in/predict/', wait = 1200):
    with open(file) as f:
        fasta = f.read()
    driver = webdriver.Chrome(driver_path,chrome_options = set_chrome_config())
    driver.implicitly_wait(wait)
    driver.get(url)
    inputElement = driver.find_element_by_xpath('//*[@id="frm1"]/p[1]/textarea')
    #xerox.copy(fasta)
    driver.execute_script("arguments[0].value = arguments[1];", inputElement,fasta)
    #inputElement.send_keys(Keys.CONTROL+ "v")
    algo = driver.find_element_by_xpath('//*[@id="frm1"]/p[6]/label/input')
    #driver.execute_script("arguments[0].click();", algo)
    algo.click()
    WebDriverWait(driver, wait).until(EC.element_to_be_clickable((By.XPATH, '//*[@id="frm1"]/p[7]/input[1]'))).click()
    #inputElement.submit()
    content = driver.page_source
    soup = bs(content,features="lxml")
    dfs = pd.read_html(soup.prettify())
    wanted = (dfs[3],dfs[4],dfs[6])
    algos = ('SVM','RFC','ANN','DAC')
    return wanted
    count = 0
    for df in wanted:
        df['Algorithm'] = algos[count]
        count+=1
    driver.close()
    final = pd.concat(wanted, ignore_index = True)
    #return  final
    
def ADAM(file, driver_path, url = 'http://bioinformatics.cs.ntou.edu.tw/ADAM/svm_tool.html'):
    with open(file) as f:
        fasta = f.read()
    driver = webdriver.Chrome(driver_path,chrome_options = set_chrome_config())
    driver.get(url)
    inputElement = driver.find_element_by_xpath('//*[@id="main2"]/form/center[1]/textarea')
    #xerox.copy(fasta)
    driver.execute_script("arguments[0].value = arguments[1];", inputElement,fasta)
    #inputElement.send_keys(Keys.CONTROL+ "v")
    inputElement.submit()
    content = driver.page_source
    soup = bs(content,features="lxml")
    driver.close()
    df = pd.read_html(soup.prettify())[1]
    df.columns = df.iloc[0]
    df.drop(df.index[0])
    return df.iloc[1: , :].reset_index(drop = True)

class ampep:

    def send(file, driver_path, mail, sleep = 2, description =  random_string(255), verbose = True):
        with open(file, 'r') as f:
            fasta = f.read()
        url = 'https://app.cbbio.online/ampep/home'
        driver = webdriver.Chrome(driver_path, options = set_chrome_config(headless = False))
        driver.get(url)
        inputElement = driver.find_element_by_xpath('//*[@id="input-65"]')
        xerox.copy(fasta)
        inputElement.send_keys(Keys.CONTROL+ "v")
        #driver.execute_script("arguments[0].value = arguments[1];", inputElement,fasta)
        driver.find_element_by_xpath('//*[@id="app"]/div/main/div/div/div/div[2]/div[2]/div/button').click()
        time.sleep(sleep)
        driver.find_element_by_xpath('//*[@id="app"]/div/main/div/div/div/div[2]/div[4]/div/button[1]/span').click()
        time.sleep(sleep)
        job_desc = driver.find_element_by_xpath('//*[@id="input-112"]')
        job_desc.send_keys(description)
        driver.find_element_by_xpath('//*[@id="app"]/div/main/div/div/div/div[2]/div[6]/div/button[1]/span').click()
        time.sleep(sleep)
        inputMail = driver.find_element_by_xpath('//*[@id="input-120"]')
        inputMail.send_keys(mail)
        time.sleep(sleep)
        submit = driver.find_element_by_xpath('//*[@id="app"]/div/main/div/div/div/div[2]/div[8]/div/button[1]/span').click()
        try:
            
            WebDriverWait(driver, 60).until(EC.presence_of_element_located((By.XPATH, '//*[@id="retrieve-page-index"]/div/div/table')))
        except TimeoutException:
            print('Timeout, took more than 60s to send job')
        if verbose == True:
            print('File sent to analysis \n to get the results call the retrieve method')
        
        driver.close()
    def retrieve(mail, driver_path, wanted = 'finished'):
        '''
        Parameters: wanted, list, default = finished
                        argument to find which job description to look. IF default = every finished job will be retrieved
        
        '''
        url = 'https://app.cbbio.online/ampep/retrieve/'+mail
        driver = webdriver.Chrome(driver_path,options = set_chrome_config())
        driver.get(url)
        time.sleep(3)
        content = driver.page_source
        soup = bs(content,features="lxml")
        driver.close()
        dfs = pd.read_html(soup.prettify())
        jobs = dfs[0]
        if wanted == 'finished':
            wanted = jobs.loc[jobs.Status == 'finished']['Job ID']
        else:
            wanted = list(jobs[jobs['Description'].isin(wanted)]['Job ID'].values)
        final = pd.DataFrame(columns = ['Unnamed: 0', 'id', 'AmPEP', 'RF-AmPEP30', 'Number of positives'])
        for ID in tqdm(range(len(wanted)), total = len(wanted)):
            ID = wanted[ID]
            url = 'https://app.cbbio.online/ampep/jobs/'+ID
            driver = webdriver.Chrome(driver_path,options = set_chrome_config())
            driver.get(url)
            time.sleep(2)
            content = driver.page_source
            soup = bs(content,features="lxml")
            df = pd.read_html(soup.prettify())[0]
            final = pd.concat([final,df],ignore_index=True)
            driver.close()
        return final