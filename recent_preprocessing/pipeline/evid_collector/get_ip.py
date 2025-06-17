import socket

#https://bluewaters.ncsa.illinois.edu/pythonnotebooks
(name, aliases, ips) = socket.gethostbyname_ex(socket.gethostname())
for ip in ips:
    print(ip)
