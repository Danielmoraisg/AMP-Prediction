from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from bs4 import BeautifulSoup as bs
import pandas as pd
import xerox
import time
from Bio import SeqIO

def portreports(file, driver_path, url = 'https://www.portoreports.com/stm'):
    with open(file) as f:
        fasta = f.read()
    driver = webdriver.Chrome(driver_path)
    driver.get(url)
    inputElement = driver.find_element_by_xpath('//*[@id="post-77"]/div/form/textarea')
    xerox.copy(fasta)
    inputElement.send_keys(Keys.CONTROL+ "v")
    inputElement.submit()
    content = driver.page_source
    soup = bs(content,features="lxml")
    driver.close()
    df = pd.read_html(soup.prettify())[0]
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
            print(f'{len(problems)*23} were not able to be used by dbaasp\nTo check sequences that were not used pleas use the argument problems = true')
        if problems:
            return df0, problems
        else:
            return df0

def campr3(file, driver_path, url = 'http://www.camp.bicnirrh.res.in/predict/'):
    with open(file) as f:
        fasta = f.read()
    driver = webdriver.Chrome(driver_path)
    driver.get(url)
    inputElement = driver.find_element_by_xpath('//*[@id="frm1"]/p[1]/textarea')
    xerox.copy(fasta)
    inputElement.send_keys(Keys.CONTROL+ "v")
    driver.find_element_by_xpath('//*[@id="frm1"]/p[6]/label/input').click()
    inputElement.submit()
    content = driver.page_source
    soup = bs(content,features="lxml")
    driver.close()
    dfs = pd.read_html(soup.prettify())
    wanted = [dfs[3],dfs[4],dfs[6]]
    algos = ['SVM','RFC','ANN','DAC']
    count = 0
    for df in wanted:
        df['Algorithm'] = algos[count]
        count+=1
    return  pd.concat(wanted, ignore_index = True)

def ADAM(file, driver_path, url = 'http://bioinformatics.cs.ntou.edu.tw/ADAM/svm_tool.html'):
    with open(file) as f:
        fasta = f.read()
    driver = webdriver.Chrome(driver_path)
    driver.get(url)
    inputElement = driver.find_element_by_xpath('//*[@id="main2"]/form/center[1]/textarea')
    xerox.copy(fasta)
    inputElement.send_keys(Keys.CONTROL+ "v")
    inputElement.submit()
    content = driver.page_source
    soup = bs(content,features="lxml")
    driver.close()
    df = pd.read_html(soup.prettify())[1]
    df.columns = df.iloc[0]
    df.drop(df.index[0])
    return df.iloc[1: , :].reset_index(drop = True)
