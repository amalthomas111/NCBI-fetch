import sys, logging

class StreamToLogger(object):
    def __init__(self, logger, level=logging.INFO):
        self.logger = logger
        self.level = level
        self.linebuf = ''

    def write(self, buf):
        #print(buf)
        for line in buf.rstrip().splitlines():
            if line!="\n" or line != "" or line !="":
                self.logger.log(self.level, line)
    def flush(self):
        pass

logging.basicConfig(
   level=logging.INFO,
   format='%(asctime)s:%(levelname)s:%(name)s:%(message)s',
   datefmt='%m-%d-%Y %I:%M:%S %p',
   filename="log.txt",
   filemode='a'
)
