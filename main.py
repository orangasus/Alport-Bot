# Don't forget to change relative paths to absolute ones, when uploading to server

# Module for working with date
import datetime

# Modules for working with Google services
import gspread
# Module for working with data
import pandas as pd
from oauth2client.service_account import ServiceAccountCredentials
from tqdm import tqdm

# My modules
import MailHandler
import PubMedHandler

# Modules for environment variables

# Setting up google API client
scope = ['https://www.googleapis.com/auth/spreadsheets', 'https://www.googleapis.com/auth/drive']
credentials = ServiceAccountCredentials.from_json_keyfile_name('google_cred.json', scope)
google_client = gspread.authorize(credentials)

# Connecting to Google sheets db
db = google_client.open("AlportBot DB Test")
sheet1 = db.worksheet('Sheet1')

# Stores all information that will be used to send emails to users
# user_report_info = [(email, formated_articles, name), ...]
user_report_info = []


def get_formated_current_date():
    today = datetime.date.today()
    res = today.strftime("%b-%d-%Y")
    month, day, year = res.split('-')
    return (day, month, year)


def get_user_info_from_row(row):
    name = 'Subscriber' if row['name'] == None else row['name']
    email = row['email']
    key_words = [x.lower() for x in row['keywords_list'].split(", ")]
    query = row['query']
    auth_limit = row['authors_limit']
    return {'name': name, 'email': email, 'key_words': key_words, 'query': query, 'auth_limit': auth_limit}


def get_prev_month_str():
    return (datetime.date.today() - datetime.timedelta(days=3)).strftime("%b-%d-%Y").split('-')[0]


if __name__ == '__main__':
    # Gets current year
    cur_year = str(datetime.datetime.now().year)
    # Gets previous month in text format e.g. Oct
    prev_month_str = get_prev_month_str()

    # counts number of sent emails just in case
    emails_sent = 0

    # Reads email credentials
    with open('email_cred.txt') as cred_txt:
        email_sender, email_password = cred_txt.readline().split(', ')

    # Reads db google sheet
    df = pd.DataFrame(sheet1.get_all_records())
    for _, row in tqdm(df.iterrows()):
        user_data = get_user_info_from_row(row)

        # fetches articles
        id_list = PubMedHandler.search(user_data['query'])['IdList']

        articles_fetch = PubMedHandler.fetch_articles_details(id_list)
        articles_sum = PubMedHandler.get_articles_summary_info(id_list)

        formated_articles = []
        month_before, month_after = None, None
        n = len(id_list)
        index = 1

        # formats fetched articles
        for i, paper in enumerate(articles_fetch['PubmedArticle']):
            if i < n - 1:
                month_before = PubMedHandler.get_formated_pubdate(articles_sum[i + 1])[1]

            d = PubMedHandler.format_article(articles_sum[i], paper, id_list[i], index, user_data['auth_limit'])

            paper_day, paper_month, paper_year = PubMedHandler.get_formated_pubdate(articles_sum[i])

            # Sometimes PubMed leaves the month field blank, this code assigns a correct value if that's the case
            if paper_month == 'Unknown' and i != 0:
                if month_after == month_before:
                    paper_month = month_after
                elif month_after != month_before:
                    if month_after != 'Unknown':
                        paper_month = month_after
                    elif month_before != 'Unknown':
                        paper_month = month_before

            # Checks if article was published in previous month
            # Then checks if users keywords are present in article's keyword list
            if paper_year == cur_year and prev_month_str == paper_month:
                try:
                    paper_keyword_list = list(paper['MedlineCitation']['KeywordList'][0])
                except:
                    paper_keyword_list = []

                paper_keyword_list = [str(el).lower() for el in paper_keyword_list]

                if any(x in user_data['key_words'] for x in paper_keyword_list):
                    formated_articles.append(d)
                    index += 1

            month_after = paper_month

        user_report_info.append((user_data['email'], formated_articles, user_data['name']))

    MailHandler.send_all_reports(user_report_info, email_sender, email_password)
