from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from bs4 import BeautifulSoup as bs
import pandas as pd
import time
from Bio import SeqIO
import os
import xerox
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
import random
import string
from tqdm import tqdm

#function to generate bathces of sequences
def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)] 

def set_chrome_config(headless = True, dest_path = 'C:/Downloads'):
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
    options.add_argument("download.default_directory="+dest_path)
    return options
    
def random_string(length):
    # With combination of lower and upper case
    result_str = ''.join(random.choice(string.ascii_letters) for i in range(length))
    # print random string
    return result_str

def portreports(file, driver_path, n_try = 100000):
    url = 'https://www.portoreports.com/stm'
    with open(file) as handle:
        records = SeqIO.parse(handle, "fasta")
        lista = list(batch(list(records), n_try))
        first = True
        for i in tqdm(range(len(lista)), total = len(lista)):
            if first:
                df0 = pd.DataFrame(columns = ['Sequence Number', 'Score', 'Prediction', 'Seq ID', 'Seq']) #starting the DataFrame
                first = False
            record = lista[i]
            SeqIO.write(record, "temp.fasta", "fasta")
            with open('temp.fasta') as f:
                fasta = f.read()
            entrys = list(SeqIO.parse('temp.fasta','fasta')) 
            os.remove("temp.fasta")
        #with open(file) as f:
        #    fasta = f.read()
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
            df['Seq ID'] = [x.id for x in entrys]
            df['Seq'] = [str(x.seq)for x in entrys]
            df0 = pd.concat([df0, df], ignore_index=True)
    return df0.reset_index(drop = True)   
        
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
            df0 = pd.DataFrame(columns = ['Seq. ID','Class','Seq ID','Seq']) #starting the DataFrame
            first = False
        problems = []
        lista = list(batch(list(records), n_try))
        for i in tqdm(range(len(lista)), total = len(lista)):
            record = lista[i]
            SeqIO.write(record, "small_temp.fasta", "fasta")
            with open('small_temp.fasta') as f:
                small_fasta = f.read()
            entrys = list(SeqIO.parse("small_temp.fasta",'fasta')) 
            os.remove("small_temp.fasta")
            driver = webdriver.Chrome(driver_path,chrome_options = set_chrome_config())
            driver.get(url)
            inputElement = driver.find_element_by_xpath('/html/body/main/div[2]/div/textarea')
            #xerox.copy(small_fasta)
            driver.execute_script("arguments[0].value = arguments[1];",inputElement,small_fasta)
            #inputElement.send_keys(Keys.CONTROL+ "v")
            driver.find_element_by_xpath('/html/body/main/div[2]/div/button').click()
            WebDriverWait(driver, 30).until(EC.visibility_of_element_located((By.XPATH, "/html/body/main/div[2]/div/table/tbody/tr/td")))
            content = driver.page_source
            soup = bs(content,features="lxml")
            driver.close()
            try:
                df = pd.read_html(soup.prettify())[0]
                df.drop(df.tail(1).index,inplace=True)
                df['Seq ID'] = [x.id for x in entrys]
                df['Seq'] = [str(x.seq)for x in entrys]
                df0 = pd.concat([df0, df], ignore_index=True)
            except IndexError:
                problems.append(record)
        print(f'{len(problems)*23} were not able to be used by dbaasp\nTo check sequences that were not used please use the argument problems = true')
    os.remove('temp.fasta')
    if problems:
        return df0, problems
    else:
        return df0
    
def campr3(file, driver_path, url = 'http://www.camp.bicnirrh.res.in/predict/', wait = 1200, n_try = 1000):
    with open(file) as handle:
        records = SeqIO.parse(handle, "fasta")
        lista = list(batch(list(records), n_try))
        first = True
        for i in tqdm(range(len(lista)), total = len(lista)):
            if first:
                df0 = pd.DataFrame(columns = ['Seq. ID.', 'Class', 'AMP Probability', 'Algorithm', 'Seq ID', 'Seq']) #starting the DataFrame
                first = False
            record = lista[i]
            SeqIO.write(record, "temp.fasta", "fasta")
            with open('temp.fasta') as f:
                fasta = f.read()
            entrys = list(SeqIO.parse('temp.fasta','fasta'))
            os.remove("temp.fasta")
    #with open(file) as f:
        #fasta = f.read()
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
            driver.close()
            dfs = pd.read_html(soup.prettify())
            wanted = (dfs[3],dfs[4],dfs[5],dfs[6])
            algos = ('SVM','RFC','ANN','DAC')
            count = 0
            #return wanted , algos , entrys, fasta
            for df in wanted:
                df['Algorithm'] = algos[count]
                try:
                    df['Seq ID'] = [x.id for x in entrys]
                    df['Seq'] = [str(x.seq)for x in entrys]
                except:
                    try:
                        used_entrys = [entrys[int(x)-1] for x in df['Seq. ID.']]
                    except:
                        return wanted , algos , entrys, fasta
                    df['Seq ID'] = [x.id for x in used_entrys]
                    df['Seq'] = [str(x.seq)for x in used_entrys]
                count+=1
            final = pd.concat(wanted, ignore_index = True)
            df0 = pd.concat([df0, final], ignore_index=True)
    return df0
    
def ADAM(file, driver_path, url = 'http://bioinformatics.cs.ntou.edu.tw/ADAM/svm_tool.html', n_try = 10000):
    with open(file) as handle:
        records = SeqIO.parse(handle, "fasta")
        lista = list(batch(list(records), n_try))
        first = True
        for i in tqdm(range(len(lista)), total = len(lista)):
            if first:
                df0 = pd.DataFrame(columns = ['Name', 'Sequence', 'Value', 'Label']) #starting the DataFrame
                first = False
            record = lista[i]
            SeqIO.write(record, "temp.fasta", "fasta")
            with open('temp.fasta') as f:
                fasta = f.read()
            entrys = list(SeqIO.parse('temp.fasta','fasta'))
            os.remove("temp.fasta")
    #with open(file) as f:
        #fasta = f.read()
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
            df = df.iloc[1: , :].reset_index(drop = True)
            df0 = pd.concat([df0, df], ignore_index=True)
    return df0.reset_index(drop = True)

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
            print(f'File {file} sent to analysis \n to get the results call the retrieve method')
        
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
def iAMP(file, driver_path, n_try = 100000):
    url = 'http://cabgrid.res.in:8080/amppred/server.php'
    with open(file) as handle:
        records = SeqIO.parse(handle, "fasta")
        lista = list(batch(list(records), n_try))
        first = True
        for i in tqdm(range(len(lista)), total = len(lista)):
            if first:
                df0 = pd.DataFrame(columns = ['name_fasta','antibacterial','antiviral','antifungal']) #starting the DataFrame
                first = False
            record = lista[i]
            SeqIO.write(record, "temp.fasta", "fasta")
            with open('temp.fasta') as f:
                fasta = f.read()
            entrys = list(SeqIO.parse('temp.fasta','fasta')) 
            os.remove("temp.fasta")
        #with open(file) as f:
        #    fasta = f.read()
            driver = webdriver.Chrome(driver_path, options = set_chrome_config())
            driver.get(url)
            inputElement = driver.find_element_by_xpath('//*[@id="form1"]/div/textarea')
            #xerox.copy(fasta)
            driver.execute_script("arguments[0].value = arguments[1];", inputElement,fasta)
            #inputElement.send_keys(Keys.CONTROL+ "v")
            inputElement.submit()
            #WebDriverWait(driver, 3).until(EC.visibility_of_element_located((By.XPATH, "/html/body/p/table/tbody/tr/td/table/tbody/tr[1]/th[2]")))
            driver.switch_to.frame(driver.find_element_by_tag_name("iframe"))
            content = driver.page_source
            soup = bs(content,features="lxml")
            driver.close()
            df = pd.read_html(soup.prettify())[0].iloc[1:]
            df.columns = df.iloc[0]
            df[['name_fasta','antibacterial','antiviral','antifungal']]
            #return df, entrys
            df = df.drop(df.index[0])
            df['Seq ID'] = [x.id for x in entrys]
            df['Seq'] = [str(x.seq)for x in entrys]
            df0 = pd.concat([df0, df], ignore_index=True)
    return df0.reset_index(drop = True)[['antibacterial', 'antiviral', 'antifungal', 'Seq ID',
       'Seq']]

def iampe(file, driver_path, n_try = 100000, download_path = 'C:/Downloads'):
    url = 'http://cbb1.ut.ac.ir/AMPClassifier/Index'
    
    with open(file) as handle:
        records = SeqIO.parse(handle, "fasta")
        lista = list(batch(list(records), n_try))
        first = True
        for i in tqdm(range(len(lista)), total = len(lista)):
            if first:
                df0 = pd.DataFrame(columns = ['Peptide sequence', ' kNN', ' SVM', ' RF', ' XGBoost','Seq ID']) #starting the DataFrame
                first = False
            record = lista[i]
            SeqIO.write(record, "temp.fasta", "fasta")
            with open('temp.fasta') as f:
                fasta = f.read()
            entrys = list(SeqIO.parse('temp.fasta','fasta')) 
            os.remove("temp.fasta")
        #with open(file) as f:
        #    fasta = f.read()
            driver = webdriver.Chrome(driver_path, options = set_chrome_config(headless = False))
            driver.get(url)
            inputElement = driver.find_element_by_xpath('//*[@id="inputSequence"]')
            #xerox.copy(fasta)
            driver.execute_script("arguments[0].value = arguments[1];", inputElement,fasta)
            #inputElement.send_keys(Keys.CONTROL+ "v")
            element = driver.find_element_by_xpath('//*[@id="classifiersForm"]/div/div[3]/div/div[2]/label[1]').click()
            element = driver.find_element_by_xpath('//*[@id="classifiersForm"]/div/div[3]/div/div[2]/label[2]').click()
            element = driver.find_element_by_xpath('//*[@id="classifiersForm"]/div/div[3]/div/div[2]/label[3]').click()
            element = driver.find_element_by_xpath('//*[@id="classifiersForm"]/div/div[3]/div/div[2]/label[4]').click()
            element = driver.find_element_by_xpath('//*[@id="classifiersForm"]/div/div[3]/div/div[2]/input[3]').click()
            element = driver.find_element_by_xpath('//*[@id="modalCl"]/div/div/div[1]/div[5]/input').click()
            file_name = driver.find_element_by_xpath('/html/body/div/div/div[2]/table/tbody/tr/td[2]/label').get_attribute("value")
            element = driver.find_element_by_xpath('/html/body/div/div/div[2]/table/tbody/tr/td[3]/input[2]').click()
            time.sleep(5)
            #WebDriverWait(driver, 3).until(EC.visibility_of_element_located((By.XPATH, "/html/body/p/table/tbody/tr/td/table/tbody/tr[1]/th[2]")))
            driver.close()
            file = download_path+r'/OutIAMPEPreds.csv'
            file = file.replace('/','//')
            df = pd.read_csv(file)
            os.remove(file)
            df['Seq ID'] = [x.id for x in entrys]
            #df['Seq'] = [str(x.seq)for x in entrys]
            df0 = pd.concat([df0, df], ignore_index=True)
    return df0.reset_index(drop = True)