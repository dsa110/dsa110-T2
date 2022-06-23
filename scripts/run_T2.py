#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# dsahead python 3.7
import T2 
import sys
import argparse
from dsautils import dsa_syslog
logger = dsa_syslog.DsaSyslogger()
logger.subsystem('software')
logger.app('T2')


def main(argv):
    parser = argparse.ArgumentParser(description='Parse input to T2 socket clients')
    parser.add_argument('--ip', type=str, default='10.41.0.46', help='ip address of heimdall', required=False)
    parser.add_argument('--ports', type=str, default='12345,12346,12347,12348', help='ports address of heimdall (comma-delimited list)', required=False)
    parser.add_argument('--trigger', type=bool, default=True, help='send trigger to dump buffer', required=False)
    parser.add_argument('--source_catalog', type=str, default='/home/ubuntu/proj/dsa110-shell/dsa110-T2/data/catalog.txt', help='set to identify triggers from sources', required=False)
    args = parser.parse_args()
    ip = args.ip
    ports = [int(port) for port in args.ports.split(',')]
    trigger = args.trigger
    source_catalog = args.source_catalog

    print(f'Running parse_socket to ip {ip} and ports {ports} with voltage trigger={trigger}')
    logger.info(f'Running parse_socket to ip {ip} and ports {ports} with voltage trigger={trigger}')
    T2.socket.parse_socket(host=ip, ports=ports, selectcols=['itime', 'idm', 'ibox', 'ibeam'], outroot="/operations/T2/cluster_output", plot_dir=None, trigger=trigger, source_catalog=source_catalog)

if __name__ == '__main__':
    main(sys.argv)
