"""Tools for interacting with the bacdive database."""
import bacdive
import logging

logger = logging.getLogger(__name__)

class BacdiveClient(bacdive.BacdiveClient):

    def getIDByNCBITaxID(self, tid):
        """Bacdive search via NCBI taxid.
        
        Parameters
        ----------
        tid : int, NCBI tax ID

        Returns
        -------
        int: number of records found
        """
        if self.public:
            baseurl = "https://bacdive.dsmz.de/"
        else:
            baseurl = "http://bacdive-dev.dsmz.local/"
        url = baseurl+"search?search=taxid:"+str(tid)
        print(url)
        resp = self.do_request(url)

        if len(resp.history) == 0:
            logger.info(f"BacdiveClient: No record found for taxid {tid}")
            return 0
        elif len(resp.history) > 0:
            raise ValueError(f"Found multiple records for taxid {tid}")
        else:
            pass

        # now get the strain ID from the response
        sid = int(resp.history[0].headers['location'].split('/')[-1])

        result = {'count': 1,
            'next': None,
            'previous': None,
            'results': [sid]}
        self.result = result
        logger.info(f"BacdiveClient: Record found for taxid {tid}")
        return 1