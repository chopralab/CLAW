<?php
// Drop-in replacement for apache_request_headers() when it's not available

if(!function_exists('apache_request_headers')) {
	function apache_request_headers() {
		static $arrHttpHeaders;
		if (!$arrHttpHeaders) {

			// Based on: http://www.iana.org/assignments/message-headers/message-headers.xml#perm-headers
			$arrCasedHeaders = array(
				// HTTP
				'Dasl'             => 'DASL',
				'Dav'              => 'DAV',
				'Etag'             => 'ETag',
				'Mime-Version'     => 'MIME-Version',
				'Slug'             => 'SLUG',
				'Te'               => 'TE',
				'Www-Authenticate' => 'WWW-Authenticate',
				// MIME
				'Content-Md5'      => 'Content-MD5',
				'Content-Id'       => 'Content-ID',
				'Content-Features' => 'Content-features',
			);
			$arrHttpHeaders = array();

			foreach($_SERVER as $strKey => $mixValue) {
				if('HTTP_' !== substr($strKey, 0, 5)) {
					continue;
				}

				$strHeaderKey = strtolower(substr($strKey, 5));

				if(0 < substr_count($strHeaderKey, '_')) {
					$arrHeaderKey = explode('_', $strHeaderKey);
					$arrHeaderKey = array_map('ucfirst', $arrHeaderKey);
					$strHeaderKey = implode('-', $arrHeaderKey);
				}
				else {
					$strHeaderKey = ucfirst($strHeaderKey);
				}

				if(array_key_exists($strHeaderKey, $arrCasedHeaders)) {
					$strHeaderKey = $arrCasedHeaders[$strHeaderKey];
				}

				$arrHttpHeaders[$strHeaderKey] = $mixValue;
			}

			/** in case you need authorization and your hosting provider has not fixed this for you:
				* VHOST-Config: 
				* FastCgiExternalServer line needs    -pass-header Authorization
				*
				* .htaccess or VHOST-config file needs: 
				  SetEnvIf Authorization "(.*)" HTTP_AUTHORIZATION=$1
				* to add the Authorization header to the environment for further processing
			 */
			if (isset($arrHttpHeaders['Authorization']) && $arrHttpHeaders['Authorization']) {
				// in case of Authorization, but the values not propagated properly, do so :)
				if (!isset($_SERVER['PHP_AUTH_USER'])) {
					list($_SERVER['PHP_AUTH_USER'], $_SERVER['PHP_AUTH_PW']) = explode(':' , base64_decode(substr($_SERVER['HTTP_AUTHORIZATION'], 6)));
				}
			}
		}
		return $arrHttpHeaders;
	}
	// execute now so other scripts have little chance to taint the information in $_SERVER
	// the data is cached, so multiple retrievals of the headers will not cause further impact on performance.
	apache_request_headers();
}