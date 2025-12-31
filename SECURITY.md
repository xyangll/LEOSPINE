# Security Policy

## Supported Versions

We provide security updates for the following versions:

| Version | Supported          |
| ------- | ------------------ |
| 1.0.x   | :white_check_mark: |
| < 1.0   | :x:                |

## Reporting a Vulnerability

We take security vulnerabilities seriously. Thank you for helping keep SPINE and our users safe.

### How to Report

**Please do not report security vulnerabilities through public GitHub issues.**

Instead, please report them via one of the following methods:

1. **GitHub Security Advisory (Recommended)**
   - Go to [https://github.com/xyangll/LEOSPINE/security/advisories/new](https://github.com/xyangll/LEOSPINE/security/advisories/new)
   - Click "Report a vulnerability"
   - Fill out the form with details about the vulnerability

2. **Private Contact**
   - Open a private security issue or contact the maintainers directly
   - Include as much detail as possible (see below)

### What to Include

When reporting a vulnerability, please provide:

- **Description**: Clear description of the vulnerability
- **Steps to Reproduce**: Detailed steps to reproduce the issue
- **Impact**: Potential impact and severity assessment
- **Environment**: 
  - Operating system and version
  - Python version
  - SPINE version
  - Branch/commit hash (if applicable)
- **Proof of Concept**: If available, include a minimal code example demonstrating the vulnerability
- **Suggested Fix**: If you have ideas for a fix, please share them

## Response Process

1. **Acknowledgment**: We will acknowledge receipt of your report within 3-5 business days
2. **Confirmation**: We will confirm and reproduce the vulnerability
3. **Assessment**: We will assess the impact and develop a fix plan
4. **Resolution**: We will work on a fix and keep you informed of progress
5. **Disclosure**: We will coordinate disclosure after a patch is released

### Timeline

- **Initial Response**: 3-5 business days
- **Status Update**: Within 2 weeks
- **Fix Timeline**: Depends on severity and complexity

## Scope

### In Scope

We are interested in security vulnerabilities in SPINE itself, including:

- Code execution vulnerabilities
- Privilege escalation
- Information disclosure
- Denial of service (DoS) affecting core functionality
- Authentication/authorization bypasses
- Input validation issues leading to security problems

### Out of Scope

The following are **not** considered security vulnerabilities:

- Known third-party CVEs (please link to the CVE instead)
- Misconfiguration risks
- Issues in dependencies (report to the dependency maintainers)
- Social engineering attacks
- Physical security issues
- Issues requiring physical access to the system

## Security Updates

Security updates will be released as patch versions (e.g., 1.0.1, 1.0.2) and will be announced in:

- GitHub Releases
- Security Advisories (if applicable)
- CHANGELOG.md

## Recognition

We appreciate responsible disclosure. With your permission, we would like to:

- Credit you in our security advisories
- Add you to our acknowledgments (if you prefer)

Thank you for helping keep SPINE secure! 🛰️

