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
    parser.add_argument('--ip', type=str, default='10.42.0.90', help='ip address of heimdall', required=False)
    parser.add_argument('--ports', type=str, default='12345,12346,12347,12348,13345,13346,13347,13348', help='ports address of heimdall (comma-delimited list)', required=False)
    #parser.add_argument('--ports', type=str, default='12345,12346,12347,12348', help='ports address of heimdall (comma-delimited list)', required=False)
    parser.add_argument('--trigger', type=bool, default=True, help='send trigger to dump buffer', required=False)
    parser.add_argument('--source_catalog', type=str, default=None, help='set to identify triggers from sources', required=False)
    #New auditing injection parameters
    parser.add_argument("--audit_injections", action="store_true", help="enable injection auditing (writes aggregated csv with injection status)")
    parser.add_argument("--audit_dump_json",action="store_true", help="In addition to aggregated csvs, dump per-injection JSON file showing live status (only used if --audit_injections is set)")
    parser.add_argument("--audit_dir", type=str, default="/operations/T2/injection_audit_results/", help="directory for dumping audit results (CSV/JSON). Ignored unless --audit_injections is set.")
    args = parser.parse_args()
    ip = args.ip
    ports = [int(port) for port in args.ports.split(',')]
    trigger = args.trigger
    source_catalog = args.source_catalog
    

    print(f'Running parse_socket to ip {ip} and ports {ports} with voltage trigger={trigger}')
    logger.info(f'Running parse_socket to ip {ip} and ports {ports} with voltage trigger={trigger}')
    T2.socket.parse_socket(host=ip, ports=ports, selectcols=['itime', 'idm', 'ibox', 'ibeam'],
                           outroot="/operations/T2/cluster_output/", plot_dir=None, trigger=trigger, source_catalog=source_catalog,
                           audit_injections=args.audit_injections, audit_dump_json=args.audit_dump_json, audit_dir=args.audit_dir)

if __name__ == '__main__':
    main(sys.argv)
