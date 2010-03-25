import smtplib
import email
import sys
import traceback
import os


# Sends an email to author if there's an exception
# Not currently enabled anywhere.
  
def send_email(text):
    
    
    SERVER = "localhost"
    
    FROM = "caton_exception_emailer@nowhere.com"
    TO = ["schulmannerism%s.com"%"@gmail"]
    
    SUBJECT = "exception by %s"%os.getlogin()
    
    TEXT = text
    
    # Prepare actual message
    
    message = """\
    From: %s
    To: %s
    Subject: %s
    
    %s
    """ % (FROM, ", ".join(TO), SUBJECT, TEXT)
    
    # Send the mail
    
    server = smtplib.SMTP(SERVER)
    server.sendmail(FROM, TO, message)
    server.quit()
    
    
def report_exception():
    traceback.print_exc(sys.exc_info())
    print "reporting exception to author John Schulman..."
    print "feel free to ctrl-C if this is taking too long"
    try:
        send_email(traceback.format_exc(sys.exc_info()))
        print "emailing succeeded"
    except Exception:
        print "emailing failed"
    
def test():
    try:
        xxx
    except Exception:
        report_exception()
    
if __name__ == '__main__':
    test()