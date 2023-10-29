from pprint import pprint

from Bio import Entrez

from ArticleScraper import AuthorScraper


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

def get_auth_per_affiliation(formated_articles):
    auth_list = formated_articles['auth_list']
    auth_per_affiliation = []
    seen_af = set()
    for auth in auth_list:
        af = auth.affiliation
        ind_af = -1
        if af in seen_af:
            for i in range(len(auth_per_affiliation)):
                if auth_per_affiliation[i][0] == af:
                    ind_af = i
                    break
            auth_per_affiliation[ind_af][1].append(auth.fullname)
        else:
            seen_af.add(af)
            auth_per_affiliation.append((af,[auth.fullname]))
    print('---------------------------->')
    pprint(auth_per_affiliation)
    return auth_per_affiliation


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


def format_article(article_sum, article_fetch, id, index, auth_limit):
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

    authScraper = AuthorScraper()
    auth_list, etAlNeeded = authScraper.scrape_article_page(link, auth_limit)
    # for el in auth_list:
    #     print(el)

    d = {'title': name, 'date': date, 'doi': doi, 'abstract': abstract, 'link': link, 'index': index,
         'auth_list': auth_list, 'et_al': etAlNeeded}
    auth_per_af = get_auth_per_affiliation(d)
    d['auth_per_af'] = auth_per_af
    return d


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
