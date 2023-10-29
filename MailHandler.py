import ssl
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from pprint import pprint

from jinja2 import Template


def format_article_for_email(article):
    result = '{}\n{}\n\n'.format(article['title'], article['date'])
    for el in article['abstract']:
        heading, text = el[0], el[1]
        result += heading + '\n' + text + '\n\n'
    result += 'DOI: {}\n{}\n\n\n\n'.format(article['doi'], article['link'])
    return result

def create_plain_text_email_content(formated_articles, name):
    num = len(formated_articles)
    intro = f'Hi, {name}! Here are some new papers for you that were published on PubMed during last month. Today I have {num} articles for you. Hope you find them interesting!\n\n\n'
    body = intro
    for article in formated_articles:
        body += format_article_for_email(article)
    return body

def create_html_email_content(formated_articles, name):
    num = len(formated_articles)
    file = open('template.html')
    tem = Template(file.read())
    msg = tem.render(formated_articles=formated_articles, num=num, name=name)
    return msg

def send_report(formated_articles, name, email, email_sender, email_password):
    num = len(formated_articles)
    subject = 'Latest Articles from PubMed'

    message = MIMEMultipart('alternative')
    message['From'] = email_sender
    message['To'] = email
    message['Subject'] = subject

    content_text = create_plain_text_email_content(formated_articles, name)
    content_html = create_html_email_content(formated_articles, name)

    plain_version = MIMEText(content_text, 'plain')
    html_version = MIMEText(content_html, 'html')
    message.attach(plain_version)
    message.attach(html_version)

    ctx = ssl.create_default_context()
    with smtplib.SMTP_SSL('smtp.gmail.com', 465, context=ctx) as smtp:
        smtp.login(email_sender, email_password)
        smtp.sendmail(email_sender, email, message.as_string())