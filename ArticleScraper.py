from pprint import pprint
from bs4 import BeautifulSoup
import requests
import lxml


class Author():
    def __init__(self, fullname, affiliation, position):
        self.fullname = str(fullname)
        self.affiliation = str(affiliation)
        self.position = int(position)

    def __str__(self):
        return f"{self.fullname}\n{self.affiliation}\n{self.position}"

    def __lt__(self, other):
        return int(self.position) < int(other.position)

class AuthorScraper():
    def scrape_article_page(self, link, limit_auth):
        response_html = requests.get(link).text
        soup = BeautifulSoup(response_html, 'lxml')

        authors_list_html = soup.find(name='div', attrs={'class': 'authors-list'})
        auth_list = authors_list_html.find_all(name='span', attrs={'class': 'authors-list-item'})

        final_auth_list = []
        for el in auth_list:
            fullname = el.find('a', attrs={'class': 'full-name'}).text
            affiliation_html = el.find('a', attrs={'class': 'affiliation-link'})
            try:
                position = affiliation_html.text.strip()
            except:
                position = 0

            try:
                affiliation = affiliation_html.get('title')
            except:
                affiliation = 'No Affiliation provided'
                print(affiliation)

            final_auth_list.append(Author(fullname, affiliation, position))
        return final_auth_list[:limit_auth], len(final_auth_list) > limit_auth


if __name__ == 'main':
    cl = AuthorScraper()
    l = cl.scrape_article_page('https://pubmed.ncbi.nlm.nih.gov/37537573/', 4)
    for el in l:
        print(el)
