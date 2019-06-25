import requests


BASE = "http://www.uniprot.org"
KB_ENDPOINT = "/uniprot/"
ID = "P12345.fasta"


if __name__ == "__main__":
    print("HI")
    fullUrl = BASE + KB_ENDPOINT + ID
    result = requests.get(fullUrl)
    if result.ok:
        print(result.text)
        fL = result.text.split("\n")
        fS = "".join(fL[1:])
        print(fS)
    else:
        print("Something went wrong ", result.status_code)
