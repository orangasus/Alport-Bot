#TODO:
# Try to make it async wherever possible
# Change email address
# Email Susie: mentioning org. name, list of emails, excel sheet?
# Upload to server until the 31st
# Upload to Git


from pprint import pprint

# Module for working with date
import datetime

# Modules for environment variables
from pathlib import Path
from dotenv import load_dotenv
import os

# Modules for working with Google services
import gspread
from oauth2client.service_account import ServiceAccountCredentials

# Module for working with data
import pandas as pd

from tqdm import tqdm

# My modules
import MailHandler
import PubMedHandler


# Setting up google API client
scope = ['https://www.googleapis.com/auth/spreadsheets', 'https://www.googleapis.com/auth/drive']
credentials = ServiceAccountCredentials.from_json_keyfile_name('google_cred.json', scope)
google_client = gspread.authorize(credentials)

# Connecting to Google sheets db
db = google_client.open("AlportBot DB")
sheet1 = db.worksheet('Sheet1')


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


def is_end_of_month():
    dt = datetime.datetime.now()
    todays_month = dt.month
    tomorrows_month = (dt + datetime.timedelta(days=1)).month
    return tomorrows_month != todays_month


if __name__ == '__main__':
    if True:
        # gets current date
        now = get_formated_current_date()

        # counts number of sent emails just in case
        emails_sent = 0

        # Loads enviroment variables from .env file
        current_dir = Path(__file__).resolve().parent if "__file__" in locals() else Path.cwd()
        envars = current_dir / ".env"
        load_dotenv(envars)

        # Reads environment variables
        email_sender = os.getenv('email_sender')
        email_password = os.getenv('email_password')

        # Reads db google sheet
        df = pd.DataFrame(sheet1.get_all_records())
        for _, row in tqdm(df.iterrows()):
            user_data = get_user_info_from_row(row)

            # fetches articles
            id_list = PubMedHandler.search(user_data['query'])['IdList']

            articles_fetch = PubMedHandler.fetch_articles_details(id_list)
            articles_sum = PubMedHandler.get_articles_summary_info(id_list)

            # formats fetched articles
            formated_articles = []
            month_before, month_after = None, None
            n = len(id_list)
            index = 1
            for i, paper in enumerate(articles_fetch['PubmedArticle']):
                if i < n - 1:
                    month_before = PubMedHandler.get_formated_pubdate(articles_sum[i + 1])[1]

                d = PubMedHandler.format_article(articles_sum[i], paper, id_list[i], index, user_data['auth_limit'])
                # print(d['auth_list'][0])
                # print('--------------------------')

                paper_day, paper_month, paper_year = PubMedHandler.get_formated_pubdate(articles_sum[i])

                if paper_month == 'Unknown' and i != 0:
                    if month_after == month_before:
                        paper_month = month_after
                    elif month_after != month_before:
                        if month_after != 'Unknown':
                            paper_month = month_after
                        elif month_before != 'Unknown':
                            paper_month = month_before

                if paper_year == now[2] and now[1] in paper_month:
                    try:
                        paper_keyword_list = list(paper['MedlineCitation']['KeywordList'][0])
                    except:
                        paper_keyword_list = []

                    paper_keyword_list = [str(el).lower() for el in paper_keyword_list]

                    if any(x in user_data['key_words'] for x in paper_keyword_list):
                        formated_articles.append(d)
                        index += 1

                month_after = paper_month

            # Sends report
            try:
                MailHandler.send_report(formated_articles, user_data['name'], user_data['email'],
                                        email_sender, email_password)
            except Exception as e:
                print(e)
            else:
                emails_sent += 1
        print(f"Sent Emails: {emails_sent}")