import smtplib
import ssl
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

from jinja2 import Template

# Formats article info for plain text email version
def format_article_for_email(article):
    result = '{}\n{}\n\n'.format(article['title'], article['date'])
    for el in article['abstract']:
        heading, text = el[0], el[1]
        result += heading + '\n' + text + '\n\n'
    result += 'DOI: {}\n{}\n\n\n\n'.format(article['doi'], article['link'])
    return result

# Creates text version of email
def create_plain_text_email_content(formated_articles, name):
    num = len(formated_articles)
    intro = f'Hi, {name}! Here are some new papers for you that were published on PubMed during last month. Today I have {num} articles for you. Hope you find them interesting!\n\n\n'
    body = intro
    for article in formated_articles:
        body += format_article_for_email(article)
    return body


# Creates html version of email
def create_html_email_content(formated_articles, name):
    num = len(formated_articles)
    file = open('template.html')
    tem = Template(file.read())
    msg = tem.render(formated_articles=formated_articles, num=num, name=name)
    return msg


# Creates a digest email message for user
def get_report(email_sender, email, formated_articles, name):
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

    return message

# Creates all emails for all users and sends them
def send_all_reports(user_report_info, email_sender, email_password):
    addressee_total = len(user_report_info)
    emails_sent = 0
    ctx = ssl.create_default_context()
    with smtplib.SMTP_SSL('smtp.gmail.com', 465, context=ctx) as smtp:
        smtp.login(email_sender, email_password)
        for el in user_report_info:
            try:
                message = get_report(email_sender, el[0], el[1], el[2])
                smtp.sendmail(email_sender, el[0], message.as_string())
                emails_sent += 1
            except:
                print(f'Unable to send email to: {el[0]}')
    print(f'Emails sent: {emails_sent} out of {addressee_total}')
