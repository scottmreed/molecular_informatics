from email.parser import Parser

# Raw email message as a string
raw_email = """\
From: sender@example.com
To: recipient@example.com
Subject: Test Email

This is the body of the email.
"""

# Parse the raw email message
email_parser = Parser()
parsed_email = email_parser.parsestr(raw_email)

# Extract details
subject = parsed_email['Subject']
from_address = parsed_email['From']
to_address = parsed_email['To']
body = parsed_email.get_payload()

# Output the parsed details
print("Subject:", subject)            # Output: Subject: Test Email
print("From:", from_address)          # Output: From: sender@example.com
print("To:", to_address)              # Output: To: recipient@example.com
print("Body:", body)                  # Output: Body: This is the body of the email.
