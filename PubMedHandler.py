from Bio import Entrez


def search(query, sort='pub_date'):
    Entrez.email = "@orangasus@gmail.com"
    handle = Entrez.esearch(db='pubmed',
                            sort=sort,
                            retmax='100',
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle)
    return results

def get_articles_summary_info(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'orangasus@gmail.com'
    handle = Entrez.esummary(db='pubmed',
                             retmode='xml',
                             id=ids)
    results = list(Entrez.read(handle))
    return results

def fetch_articles_details(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'orangasus@gmail.com'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results

def format_article_abstract(abstract):
    formated_abstract = []
    for line in abstract:
        content = str(line)
        try:
            heading = line.attributes['Label']
        except:
            heading = ''
        formated_abstract.append((heading, content))
    return formated_abstract

def format_article(article_sum, article_fetch, id, index):
    try:
        doi = article_sum['DOI']
    except:
        doi = 'Unknown'
    name = article_sum['Title']
    date = article_sum['PubDate']
    link = f'https://pubmed.ncbi.nlm.nih.gov/{id}/'
    try:
        abstract = format_article_abstract(article_fetch['MedlineCitation']['Article']['Abstract']['AbstractText'])
    except:
        abstract = []

    return {'title': name, 'date': date, 'doi': doi, 'abstract': abstract, 'link': link, 'index' : index}

def get_formated_pubdate(article_sum):
    initial = article_sum['PubDate']
    paper_year, paper_month, paper_day = None, None, None
    try:
        l_date = initial.split(' ')
        paper_year = l_date[0]
        paper_month = l_date[1]
        paper_day = l_date[2]
    except:
        if paper_month == None:
            paper_month = 'Unknown'
        if paper_day == None:
            paper_day = 'Unknown'

    return (paper_day, paper_month, paper_year)