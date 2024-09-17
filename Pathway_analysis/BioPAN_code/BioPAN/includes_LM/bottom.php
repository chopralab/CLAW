</div>
		</div>

		<footer class="bg-dark" style="padding-top:1rem;">
			<div class="footer-bottom bg-light text-muted" style='padding-top:0.4rem; padding-bottom:0rem; text-align: center;'>
				<div class="container">
					<div class="row">
						<div class="col" style='text-align: center;'>
							<p style='font-size:1.6rem;'>Participating Institutions</p>
						</div>
					</div>
				</div>
			</div>

			<div class='bg-light' style='padding-top:0.4rem;'>
				<div class="container">
					<div class="row align-items-end" style='border-bottom-style:solid; border-color: #cccccc; padding-bottom:2rem;'>
						<div class="col-sm" style='text-align:center; padding:1rem;'>
							<img src="/includes_LM/images/cardiff_logo.png" style="width: 4rem;" alt="Cardiff University logo">
							<br><a href='https://www.cardiff.ac.uk'>Cardiff University</a>
						</div>
						<div class="col-sm" style='text-align:center; padding:1rem;'>
							<img src="/includes_LM/images/BI-2016-small.png" style="width: 4.9rem;" alt="Babraham Institutei logo">
							<br><a href='https://www.babraham.ac.uk/'>Babraham Institute</a>
						</div>
						<div class="col-sm" style='text-align:center; padding:1rem;'>
							<img src="/includes_LM/images/UCSanDiegoLogo-BlueGold.png" style="width: 150px;" alt="University of California, San Diego logo">
							<br><a href='https://www.ucsd.edu'>University of California, San Diego</a>
						</div>
					</div>
				</div>
			</div>

			<div class="footer-bottom bg-light text-muted" style='margin-top:0.0rem; padding-top:1.4rem; padding-bottom:0rem; text-align: center;'>
				<div class="container">
					<div class="row">
						<div class="col" style='text-align: center;'>
							<p style='font-size:1.6rem;'>Funded by</p>
						</div>
					</div>
				</div>
			</div>

			<div class="footer-bottom bg-light text-muted" style='padding-top:0rem; text-align: center;'>
				<div class="container">
					<div class="row">
						<div class="col" style='text-align: center;'>
							<p>
								<img src="/includes_LM/images/wellcome-logo-black.jpg" style="width: 4rem;" alt="Wellcome Trust logo">
								<br>
								<a href='https://wellcome.ac.uk/'>Wellcome Trust</a>
							</p>
						</div>
					</div>
				</div>
			</div>

			<div class="m-t-3">
			</div>

			<div class="footer-top p-y-2">
				<div class="container-fluid">
				</div>
			</div>

			<div class="footer-bottom text-muted" style='background-color:black; padding-top:0.4rem; padding-right:0rem;'>
				<div class="container" style='width:100%;'>
					<div class="row" style='width:100%;'>

						<div class="col-sm-4">
							<span class='text-muted small'>
								&copy; 2003-<?= (int)date('Y') ?> LIPID MAPS® Lipidomics Gateway
							</span>
						</div>

						<div class="col-sm-2">
							<span class='text-muted small'>
								<a href="http://www.w3.org/html/logo/" class="no-icon">
									<img src="https://www.w3.org/html/logo/badge/html5-badge-h-css3-semantics.png" height='30' alt="HTML5 Powered with CSS3 / Styling, and Semantics" title="HTML5 Powered with CSS3 / Styling, and Semantics">
								</a>
							</span>
						</div>


						<div class="col-sm">
							<ul class="list-inline text-muted small">
								<li class="list-inline-item"><a href="/about/termsofuse.php">Terms of Use</a></li>
								<li class="list-inline-item"><a href="/about/howtocite.php">How to Cite</a></li>
								<li class="list-inline-item"><a href="/about/howtolink.php">How to Link</a></li>
								<li class="list-inline-item"><a href="/about/browsers.php">Browsers Supported</a></li>
								<li class="list-inline-item"><a href="/about/cookies.php">Manage Cookies</a></li>
							</ul>
						</div>
					</div>
				</div>
			</div>
		</footer>

		<!-- Cookie Notice Modal -->
		<div class="modal fade" id="cookie_notice_modal" tabindex="-1" role="dialog" aria-labelledby="cookie_notice_label" aria-hidden="true">
			<div class="modal-dialog modal-dialog-centered modal-lg" role="document">
				<div class="modal-content">
					<div class="modal-body">
						<div class="modal-title h2" id="cookie_notice_label">We use Cookies</div>
						<p class="my-3">
							This site uses cookies and other tracking technologies to assist with navigation and your ability to provide feedback, analyze your use of our products and services, assist with our promotional and marketing efforts, and provide content from third parties.
						</p>
						<div class="text-right">
							<button id="cookie_notice_prefs" type="button" class="btn btn-link">More Options</button>
							<button id="cookie_notice_accept" type="button" class="btn btn-primary">Accept</button>
						</div>
					</div>
				</div>
			</div>
		</div>

		<script type="text/javascript" src="js/librairies/jquery.min.js"></script>
		<script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.0/umd/popper.min.js" integrity="sha384-cs/chFZiN24E4KMATLdqdvsezGxaGsi4hLGOzlXwp5UZB1LY//20VyM2taTB4QvJ" crossorigin="anonymous"></script>
		<script src="https://stackpath.bootstrapcdn.com/bootstrap/4.1.0/js/bootstrap.min.js" integrity="sha384-uefMccjFJAIv6A+rW+L4AHf99KvxDjWSu1z9VI8SKNVmz4sk7buKt/6v9KI65qnm" crossorigin="anonymous"></script>
		<script src="https://cdn.jsdelivr.net/npm/js-cookie@2/src/js.cookie.min.js"></script>
	</body>
</html>