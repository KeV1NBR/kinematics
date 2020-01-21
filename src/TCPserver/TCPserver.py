import socket, sys
from thread import *

class Server:
    def __init__(self, clientTask):
        try:
            self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.clientTask = clientTask
        except socket.error, msg:
            sys.stderr.write("[ERROR] %s\n" % msg[1])
            sys.exit(1)


    def setting(self, portName, connetLimit):
        self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.sock.bind(('', portName))
        self.sock.listen(connetLimit)

    def open(self):
        while True:
            (csock, adr) = self.sock.accept()
            print "Client Info: ", csock, adr
            start_new_thread(self.clientTask, (csock,))

    def close(self):
        self.sock.close()

